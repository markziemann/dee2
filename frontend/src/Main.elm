module Main exposing (..)

import Browser exposing (Document)
import Browser.Navigation as Nav
import Elastic exposing (parse, serialize)
import Helpers exposing (decodeSearchResults)
import HomePage.Main as HPMain
import Html exposing (..)
import Html.Attributes exposing (..)
import Http
import Info exposing (introduction)
import Nav exposing (navBar)
import ResultsPage.Main as RPMain exposing (newSearchResults)
import ResultsPage.Types exposing (MaybeExpired(..))
import ResultsPage.Views exposing (viewSearchResults)
import Routes
import SearchPage.Helpers exposing (getMode, updateLevel)
import SearchPage.Main as SPMain
import SearchPage.Types exposing (Level(..), Mode(..), SearchParameters(..))
import SearchPage.Views
import SharedTypes exposing (PaginationOffset, RemoteData(..), toWebData)
import Types exposing (..)
import Url
import Url.Builder


init : () -> Url.Url -> Nav.Key -> ( Model, Cmd Msg )
init flags url navKey =
    let
        route =
            Routes.determinePage url

        model =
            { navKey = navKey
            , url = url
            , homePage = HPMain.init navKey
            , searchPage = SPMain.init navKey
            , resultsPage = RPMain.init navKey Nothing (Current NotAsked)
            , route = route
            }
    in
    case model.route of
        Routes.ResultsRoute searchUrlParameters ->
            ( model, search searchUrlParameters )

        _ ->
            ( model, Cmd.none )



---- UPDATE ----


search : SearchParameters -> Cmd Msg
search ((SearchParameters searchLevel searchMode searchString paginationOffset) as searchParameters) =
    let
        endPoint =
            "api/search/"

        parsedSearchString =
            case searchMode of
                Strict ->
                    Result.withDefault (Elastic.Word searchString) (parse searchString)
                        |> serialize

                Fuzzy ->
                    Elastic.Word searchString
                        |> serialize
    in
    Http.get
        { url = endPoint ++ (Url.Builder.toQuery <| Routes.searchResultParams searchParameters)
        , expect =
            Http.expectJson
                (GotHttpSearchResponse searchParameters << toWebData)
                (decodeSearchResults paginationOffset)
        }


update : Msg -> Model -> ( Model, Cmd Msg )
update msg model =
    let
        fromHomePage =
            \( mdl, cmd ) ->
                ( { model | homePage = mdl }, Cmd.map GotHomePageMsg cmd )

        fromSearchPage =
            \( mdl, cmd ) ->
                ( { model | searchPage = mdl }, Cmd.map GotSearchRunsPageMsg cmd )

        fromResultsPage =
            \( mdl, cmd ) ->
                ( { model | resultsPage = mdl }, Cmd.map GotResultsPageMsg cmd )
    in
    case msg of
        GotHomePageMsg message ->
            HPMain.update message model.homePage
                |> fromHomePage

        GotSearchRunsPageMsg message ->
            SPMain.update message model.searchPage
                |> fromSearchPage

        GotSearchProjectsPageMsg message ->
            ( model, Cmd.none )

        GotResultsPageMsg message ->
            RPMain.update message model.resultsPage
                |> fromResultsPage

        LinkClicked urlRequest ->
            case urlRequest of
                Browser.Internal url ->
                    ( model, Nav.pushUrl model.navKey (Url.toString url) )

                Browser.External href ->
                    ( model, Nav.load href )

        GotHttpSearchResponse searchParameters webDataSearchResults ->
            ( { model | resultsPage = newSearchResults model.resultsPage searchParameters webDataSearchResults }
            , Cmd.none
            )

        UrlChanged url ->
            let
                route =
                    Routes.determinePage url

                newModel =
                    { model | url = url, route = route }
            in
            case route of
                -- Any page that is associated with a url should be handled in this case statement.
                -- While it is possible to handle this logic else where in the code base
                -- doing so will be trouble for handling browser forward and backward
                -- events. Therefor any action that needs to change the page should simply
                -- push a new url which will be caught and handled here!
                Routes.ResultsRoute searchParameters ->
                    ( { newModel
                        | resultsPage = newSearchResults model.resultsPage searchParameters Loading
                      }
                    , search searchParameters
                    )

                Routes.SearchRunsRoute ->
                    ( { newModel | searchPage = updateLevel Runs model.searchPage }, Cmd.none )

                Routes.SearchProjectsRoute ->
                    ( { newModel | searchPage = updateLevel Projects model.searchPage }, Cmd.none )

                _ ->
                    ( newModel, Cmd.none )



---- VIEW ----


pageLayout : List (Html Msg) -> List (Html Msg)
pageLayout content =
    [ Html.map GotHomePageMsg navBar
    , div [ class "container my-5 mx-auto" ] content
    , introduction
    ]


pageView : Model -> List (Html Msg)
pageView model =
    let
        fromHomePage =
            List.map (Html.map GotHomePageMsg)

        fromSearchPage =
            List.map (Html.map GotSearchRunsPageMsg)

        fromResultsPage =
            List.map (Html.map GotResultsPageMsg)
    in
    pageLayout <|
        case model.route of
            Routes.HomeRoute ->
                fromHomePage [ HPMain.view ]

            Routes.SearchRunsRoute ->
                fromSearchPage
                    [ SearchPage.Views.view model.searchPage ]

            Routes.SearchProjectsRoute ->
                fromSearchPage [ SearchPage.Views.view model.searchPage ]

            Routes.ResultsRoute (SearchParameters _ _ _ paginationOffset) ->
                fromResultsPage
                    (viewSearchResults model.resultsPage paginationOffset)

            Routes.Unknown ->
                [ text "Hmm... I don't recognise that url." ]


view : Model -> Document Msg
view model =
    { title = "Digital Expression Explorer 2"
    , body = pageView model
    }


subscriptions : Model -> Sub Msg
subscriptions model =
    case model.route of
        Routes.HomeRoute ->
            Sub.none

        Routes.ResultsRoute _ ->
            Sub.none

        Routes.Unknown ->
            subscriptions { model | route = Routes.HomeRoute }

        Routes.SearchRunsRoute ->
            Sub.batch
                [ Sub.map GotSearchRunsPageMsg <| SPMain.subscriptions model.searchPage
                ]

        Routes.SearchProjectsRoute ->
            Sub.none



---- PROGRAM ----


main : Program () Model Msg
main =
    Browser.application
        { view = view
        , init = init
        , update = update
        , subscriptions = subscriptions
        , onUrlRequest = LinkClicked
        , onUrlChange = UrlChanged
        }
