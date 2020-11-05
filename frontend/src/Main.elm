port module Main exposing (..)

import Array
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
import Table
import Types exposing (..)
import Url
import Url.Builder as UBuilder


port consoleLog : String -> Cmd msg


init : () -> Url.Url -> Nav.Key -> ( Model, Cmd Msg )
init flags url navKey =
    ( { navKey = navKey
      , url = url
      , searchPage = SPMain.init
      , resultsPage = RPMain.init
      , page = Routes.determinePage url
      }
    , Cmd.none
    )



---- UPDATE ----


updateUrl : Nav.Key -> Cmd msg -> SearchPage.Types.SearchOutMsg -> Cmd msg
updateUrl navKey cmd outMsg =
    Cmd.batch
        [ cmd
        , UBuilder.absolute [ Routes.searchResultsSlug ] [ UBuilder.string "q" outMsg.searchString ]
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
            \( mdl, cmd, maybeSearchOutMsg ) ->
                ( { model
                    | searchPage = mdl
                    , resultsPage = maybeUpdate updateSearchData maybeSearchOutMsg model.resultsPage
                  }
                , Cmd.map GotSearchPageMsg cmd |> maybeAddCommand (updateUrl model.navKey) maybeSearchOutMsg
                )

        fromResultsPage =
            \( mdl, cmd, maybeResultsOutMsg ) ->
                ( { model | resultsPage = mdl }
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

        LinkClicked urlRequest ->
            case urlRequest of
                Browser.Internal url ->
                    ( model, Nav.pushUrl model.navKey (Url.toString url) )

                Browser.External href ->
                    ( model, Nav.load href )

        UrlChanged url ->
            ( { model | url = url, page = Routes.determinePage url }
            , Cmd.none
            )



---- VIEW ----


pageLayout : List (Html Msg) -> List (Html Msg)
pageLayout content =
    [ navbar
    , div [ class "d-flex justify-content-center" ] [ div [ class "container my-5 no-gutters" ] content ]
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
        case model.page of
            Routes.HomePage pageData ->
                fromSearchPage
                    [ viewLargeSearchBar model.searchPage
                    , viewSearchModeSelector model.searchPage.searchMode
                    , viewSearchButton
                    ]

            Routes.SearchResultsPage pageData ->
                fromResultsPage
                    (viewSearchResults model.resultsPage)


view : Model -> Document Msg
view model =
    { title = "Digital Expression Explorer 2"
    , body = pageView model
    }


subscriptions : Model -> Sub Msg
subscriptions model =
    case model.page of
        Routes.HomePage route ->
            Sub.batch
                [ Sub.map GotSearchPageMsg <| SPMain.subscriptions model.searchPage
                ]

        Routes.SearchResultsPage route ->
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
