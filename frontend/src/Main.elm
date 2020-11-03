port module Main exposing (..)

import Array
import Browser exposing (Document)
import Browser.Navigation as Nav
import Html exposing (..)
import Html.Attributes exposing (..)
import Info exposing (introduction)
import MainTypes exposing (..)
import MainViews exposing (viewSearchResults)
import Nav exposing (navbar)
import Routes
import SearchBar
import SearchBarHelpers exposing (delay)
import SearchBarTypes
import SearchBarViews exposing (..)
import Table
import Url


port consoleLog : String -> Cmd msg


init : () -> Url.Url -> Nav.Key -> ( Model, Cmd Msg )
init flags url navKey =
    ( { navKey = navKey
      , url = url
      , searchBar = SearchBar.init
      , page = Routes.determinePage url
      , searchHits = Nothing
      , searchResultRows = Nothing
      , resultsTableState = Table.initialSort "id"
      , resultsTableQuery = ""
      , downloading = False
      }
    , Cmd.none
    )



---- UPDATE ----


updateSearchData : Model -> SearchBarTypes.OutMsg -> Model
updateSearchData model outMsg =
    { model
        | searchHits = Just outMsg.hits
        , searchResultRows = Just outMsg.rows
    }


updateUrl : Nav.Key -> Cmd msg -> SearchBarTypes.OutMsg -> Cmd msg
updateUrl navKey cmd outMsg =
    Cmd.batch
        [ cmd
        , [ Routes.searchResultsSlug, outMsg.searchString ]
            |> String.join "/?q="
            |> Nav.pushUrl navKey
        ]


maybeUpdate : (Model -> a -> Model) -> Maybe a -> Model -> Model
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
            \( mdl, cmd, searchResults ) ->
                ( { model | searchBar = mdl } |> maybeUpdate updateSearchData searchResults
                , Cmd.map GotSearchBarMsg cmd |> maybeAddCommand (updateUrl model.navKey) searchResults
                )

        noOp =
            ( model, Cmd.none )
    in
    case msg of
        GotSearchBarMsg message ->
            SearchBar.update message model.searchBar
                |> fromSearchBar

        ResultClicked result ->
            let
                searchResultRows =
                    Maybe.map (Array.set result.id { result | selected = not result.selected })
                        model.searchResultRows
            in
            ( { model | searchResultRows = searchResultRows }, Cmd.none )

        SetResultsTableQuery resultsTableQuery ->
            ( { model | resultsTableQuery = resultsTableQuery }
            , Cmd.none
            )

        SetResultsTableState resultsTableState ->
            ( { model | resultsTableState = resultsTableState }
            , Cmd.none
            )

        DownloadRequested ->
            ( { model | downloading = True }, delay 5000 DownloadButtonReset )

        DownloadButtonReset ->
            ( { model | downloading = False }, Cmd.none )

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
    , div [ class "container my-5 mx-5 mx-auto" ] content
    , introduction
    ]


pageView : Model -> List (Html Msg)
pageView model =
    let
        fromSearchBar =
            Html.map GotSearchBarMsg
    in
    pageLayout <|
        case model.page of
            Routes.HomePage pageData ->
                List.map fromSearchBar
                    [ viewLargeSearchBar model.searchBar
                    , viewSearchModeSelector model.searchBar.searchMode
                    , viewSearchButton
                    ]

            Routes.SearchResultsPage pageData ->
                viewSearchResults model


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
                [ Sub.map GotSearchBarMsg <| SearchBar.subscriptions model.searchBar
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
