module Main exposing (..)

import Browser exposing (Document)
import Browser.Navigation as Nav
import Html exposing (..)
import Html.Attributes exposing (..)
import Info exposing (introduction)
import Nav exposing (navbar)
import ResultsPage.Helpers exposing (updateSearchData)
import ResultsPage.Main as RPMain
import ResultsPage.Views exposing (viewSearchResults)
import Routes
import SearchPage.Main as SPMain
import SearchPage.Types
import SearchPage.Views exposing (viewLargeSearchBar, viewSearchButton, viewSearchModeSelector)
import SharedTypes exposing (PaginationOffset)
import Types exposing (..)
import Url


init : () -> Url.Url -> Nav.Key -> ( Model, Cmd Msg )
init flags url navKey =
    let

        route = Routes.determinePage url

        model =
            { navKey = navKey
            , url = url
            , searchPage = SPMain.init
            , resultsPage = RPMain.init
            , route = route
            , defaultPaginationOffset = PaginationOffset 20 0
            }
    in
    case model.route of
        Routes.SearchRoute { searchMode, searchString, paginationOffset } ->
            update (RequestSearch searchMode searchString paginationOffset) model

        _ ->
            ( model, Cmd.none )


resetModel : Model -> Model
resetModel model =
    { model
        | searchPage = SPMain.init
        , resultsPage = RPMain.init
    }



---- UPDATE ----


updateUrl : Nav.Key -> Cmd msg -> SearchPage.Types.OutMsg -> Cmd msg
updateUrl navKey cmd outMsg =
    Cmd.batch
        [ cmd
        , Routes.searchResultsRoute outMsg.searchMode outMsg.searchString outMsg.paginationOffset
            |> Nav.pushUrl navKey
        ]


maybeUpdate : (model -> a -> model) -> Maybe a -> model -> model
maybeUpdate updateFunc maybe model =
    case maybe of
        Just value ->
            updateFunc model value

        Nothing ->
            model


maybeAddCommand : (Cmd msg -> a -> Cmd msg) -> Maybe a -> Cmd msg -> Cmd msg
maybeAddCommand addCommandFunc maybe cmd =
    case maybe of
        Just value ->
            addCommandFunc cmd value

        Nothing ->
            cmd


update : Msg -> Model -> ( Model, Cmd Msg )
update msg model =
    let
        fromSearchBar =
            \( mdl, cmd, maybeOutMsg ) ->
                ( { model
                    | searchPage = mdl
                    , resultsPage = maybeUpdate updateSearchData maybeOutMsg model.resultsPage
                  }
                , Cmd.map GotSearchPageMsg cmd |> maybeAddCommand (updateUrl model.navKey) maybeOutMsg
                )

        fromResultsPage =
            \( mdl, cmd, maybeOutMsg ) ->
                case maybeOutMsg of
                    Just pagination ->
                        update
                            (RequestSearch
                                model.searchPage.searchMode
                                model.searchPage.searchString
                                pagination
                            )
                            { model | resultsPage = mdl }
                            |> (\( m, c ) -> ( m, Cmd.batch [ c, Cmd.map GotResultsPageMsg cmd ] ))

                    Nothing ->
                        ( { model
                            | resultsPage = mdl
                          }
                        , Cmd.map GotResultsPageMsg cmd
                        )
    in
    case msg of
        GotSearchPageMsg message ->
            SPMain.update message model.searchPage
                |> fromSearchBar

        GotResultsPageMsg message ->
            RPMain.update message model.resultsPage
                |> fromResultsPage

        RequestSearch searchMode searchString pagination ->
            SPMain.update
                (SPMain.search searchMode searchString pagination)
                model.searchPage
                |> fromSearchBar

        LinkClicked urlRequest ->
            case urlRequest of
                Browser.Internal url ->
                    ( model, Nav.pushUrl model.navKey (Url.toString url) )

                Browser.External href ->
                    ( model, Nav.load href )

        UrlChanged url ->
            let
                newModel =
                    { model | url = url, route = Routes.determinePage url }
            in
            case Routes.determinePage url of
                Routes.HomeRoute ->
                    ( resetModel newModel, Cmd.none )

                _ ->
                    ( newModel, Cmd.none )



---- VIEW ----


pageLayout : List (Html Msg) -> List (Html Msg)
pageLayout content =
    [ navbar
    , div [ class "container my-5 mx-auto" ] content
    , introduction
    ]


pageView : Model -> List (Html Msg)
pageView model =
    let
        fromSearchPage =
            List.map (Html.map GotSearchPageMsg)

        fromResultsPage =
            List.map (Html.map GotResultsPageMsg)
    in
    pageLayout <|
        case model.route of
            Routes.HomeRoute ->
                fromSearchPage
                    [ viewLargeSearchBar model.searchPage
                    , viewSearchModeSelector model.searchPage.searchMode
                    , viewSearchButton
                        model.searchPage
                        model.defaultPaginationOffset
                    ]

            Routes.SearchRoute { searchString, paginationOffset } ->
                fromResultsPage
                    (viewSearchResults model.resultsPage)

            Routes.Unknown ->
                [text "Hmm... I don't recognise that url."]


view : Model -> Document Msg
view model =
    { title = "Digital Expression Explorer 2"
    , body = pageView model
    }


subscriptions : Model -> Sub Msg
subscriptions model =
    case model.route of
        Routes.HomeRoute ->
            Sub.batch
                [ Sub.map GotSearchPageMsg <| SPMain.subscriptions model.searchPage
                ]

        Routes.SearchRoute _ ->
            Sub.none

        Routes.Unknown ->
            subscriptions { model | route = Routes.HomeRoute }



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
